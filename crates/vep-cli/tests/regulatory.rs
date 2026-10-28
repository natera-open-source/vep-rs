// Copyright (c) 2026-present Natera, Inc.
// Licensed under the Apache License, Version 2.0 (see LICENSE).

use serde_json::Value;
use std::fs;
use std::path::Path;
mod common;

fn cache(dir: &Path) -> std::path::PathBuf {
    common::regulatory_cache(dir)
}

fn annotate(dir: &Path, cache: &Path, format: &str, regulatory: bool) -> String {
    annotate_with_args(dir, cache, format, regulatory, &[])
}

fn annotate_with_args(
    dir: &Path,
    cache: &Path,
    format: &str,
    regulatory: bool,
    extra_args: &[&str],
) -> String {
    let output = dir.join(format!("out-{format}-{regulatory}"));
    let mut command = common::isolated_vep(dir);
    command
        .args([
            "--json_cache",
            cache.to_str().unwrap(),
            "--assembly",
            "GRCh38",
            "--no_stats",
            "--no_fasta",
            "--force_overwrite",
            "--format",
            "vcf",
            "--output_format",
            if format == "default" { "vep" } else { format },
        ])
        .arg("--input_file")
        .arg(dir.join("in.vcf"))
        .arg("--output_file")
        .arg(&output);
    if regulatory {
        command.arg("--regulatory");
    }
    command.args(extra_args);
    let result = command.output().unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(!String::from_utf8_lossy(&result.stderr).contains("no regulatory data"));
    fs::read_to_string(output).unwrap()
}

#[test]
fn regulatory_rows_and_extras_reach_every_text_format() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache(dir.path());
    fs::write(dir.path().join("in.vcf"), "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t106\tforward\tC\tG\t.\t.\t.\n21\t115\treverse\tA\tG\t.\t.\t.\n").unwrap();
    let json = annotate(dir.path(), &cache, "json", true);
    let records: Vec<Value> = json
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();
    let first = &records[0];
    assert_eq!(
        first["regulatory_feature_consequences"][0]["regulatory_feature_id"],
        "ENSR_FIXTURE"
    );
    assert_eq!(
        first["regulatory_feature_consequences"][0]["biotype"],
        "enhancer"
    );
    assert!(first["regulatory_feature_consequences"][0]
        .get("strand")
        .is_none());
    let motif = &first["motif_feature_consequences"][0];
    assert_eq!(motif["motif_feature_id"], "ENSM_FORWARD");
    assert_eq!(motif["motif_name"], "ENSPFM_FORWARD");
    assert_eq!(motif["motif_pos"], 2);
    assert_eq!(
        motif["transcription_factors"],
        serde_json::json!(["TF_A", "TF_B"])
    );
    assert!(motif.get("motif_score_change").is_none());
    assert!(motif.get("high_inf_pos").is_none());
    assert_eq!(first["most_severe_consequence"], "TF_binding_site_variant");
    assert_eq!(records[1]["motif_feature_consequences"][0]["motif_pos"], 4);
    assert_eq!(
        records[1]["motif_feature_consequences"][0]["transcription_factors"],
        serde_json::json!([])
    );
    for format in ["default", "tab", "vcf"] {
        let text = annotate(dir.path(), &cache, format, true);
        let rows = text
            .lines()
            .filter(|line| !line.starts_with('#'))
            .collect::<Vec<_>>()
            .join("\n");
        assert!(
            rows.contains("regulatory_region_variant"),
            "{format}: {rows}"
        );
        assert!(rows.contains("TF_binding_site_variant"), "{format}: {rows}");
        assert!(rows.contains("enhancer"), "{format}: {rows}");
        assert!(rows.contains("ENSPFM_FORWARD"), "{format}: {rows}");
        assert!(rows.contains("TF_A"), "{format}: {rows}");
    }
    let disabled = annotate(dir.path(), &cache, "json", false);
    for line in disabled.lines() {
        let record: Value = serde_json::from_str(line).unwrap();
        assert!(record.get("regulatory_feature_consequences").is_none());
        assert!(record.get("motif_feature_consequences").is_none());
        assert_eq!(record["most_severe_consequence"], "intergenic_variant");
    }
    // Metadata is authoritative even when a regulatory tree exists.
    fs::write(cache.join("info.json"), r#"{"assembly":"GRCh38"}"#).unwrap();
    fs::write(cache.join("regulatory/21/1-1000000.json"), "[").unwrap();
    let mut command = common::isolated_vep(dir.path());
    // An output of STDOUT writes its warnings file beside the working directory.
    let result = command
        .current_dir(dir.path())
        .args([
            "--json_cache",
            cache.to_str().unwrap(),
            "--regulatory",
            "--no_fasta",
            "--no_stats",
            "--json",
            "--input_file",
        ])
        .arg(dir.path().join("in.vcf"))
        .args(["--output_file", "STDOUT"])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert!(String::from_utf8_lossy(&result.stderr).contains("this cache holds no regulatory data"));
    assert!(!String::from_utf8_lossy(&result.stdout).contains("motif_feature_consequences"));
}

/// Two records of the gnomAD v4.1 chromosome 21 dataset against the slice of
/// the converted release 116 cache (`fixtures/regulatory/chr21-slice.json`):
/// `rs868224642` (`21:46699255 G>C`) inside the enhancer `ENSR21_C6M55` and
/// its four `NRF1` motifs on both strands, and `21:5090724 A>T` inside two
/// motifs of different matrices and no regulatory feature. Every row is held
/// whole, in the default, tab, VCF and JSON formats, to the row Ensembl VEP
/// 116.2 writes for the same record on the same cache, with the two fields
/// vep-rs does not write, `HIGH_INF_POS` and `MOTIF_SCORE_CHANGE`, taken out
/// (`Y`, `-0.094`, `Y`, `-0.094`, `N`, `-0.070`, `N`, `-0.082` on the NRF1
/// rows; `Y`, `-0.089` and `N`, `0.039` on the two others).
#[test]
fn regulatory_rows_are_byte_exact_in_every_format() {
    let dir = tempfile::tempdir().unwrap();
    let cache = common::regulatory_cache_of_real_features(dir.path());
    fs::write(
        dir.path().join("in.vcf"),
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t5090724\t.\tA\tT\t.\t.\t.\n\
         21\t46699255\trs868224642\tG\tC\t.\t.\t.\n",
    )
    .unwrap();
    let rows = |format: &str| -> Vec<String> {
        annotate(dir.path(), &cache, format, true)
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(str::to_owned)
            .collect()
    };
    let nrf1 = |id: &str, strand: &str, pos: &str| {
        format!(
            "rs868224642\t21:46699255\tC\t-\t{id}\tMotifFeature\tTF_binding_site_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND={strand};MOTIF_NAME=ENSPFM0427;MOTIF_POS={pos};TRANSCRIPTION_FACTORS=NRF1"
        )
    };
    assert_eq!(
        rows("default"),
        [
            "21_5090724_A/T\t21:5090724\tT\t-\tENSM00000574814\tMotifFeature\tTF_binding_site_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;MOTIF_NAME=ENSPFM0290;MOTIF_POS=7;TRANSCRIPTION_FACTORS=HNF4A,NR2F1,RXRA".to_string(),
            "21_5090724_A/T\t21:5090724\tT\t-\tENSM00000574815\tMotifFeature\tTF_binding_site_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=-1;MOTIF_NAME=ENSPFM0484;MOTIF_POS=15;TRANSCRIPTION_FACTORS=RUNX3".to_string(),
            "21_5090724_A/T\t21:5090724\tT\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER".to_string(),
            "rs868224642\t21:46699255\tC\t-\tENSR21_C6M55\tRegulatoryFeature\tregulatory_region_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;BIOTYPE=enhancer".to_string(),
            nrf1("ENSM00000588223", "-1", "5"),
            nrf1("ENSM00000588224", "1", "8"),
            nrf1("ENSM00000588225", "-1", "11"),
            nrf1("ENSM00000588226", "1", "2"),
            "rs868224642\t21:46699255\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER".to_string(),
        ]
    );
    let tab_motif = |uploaded: &str,
                     location: &str,
                     allele: &str,
                     id: &str,
                     strand: &str,
                     matrix: &str,
                     pos: &str,
                     factors: &str| {
        format!(
            "{uploaded}\t{location}\t{allele}\t-\t{id}\tMotifFeature\tTF_binding_site_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t{strand}\t-\t-\t{matrix}\t{pos}\t-\t-\t{factors}"
        )
    };
    assert_eq!(
        rows("tab"),
        [
            tab_motif("21_5090724_A/T", "21:5090724", "T", "ENSM00000574814", "1", "ENSPFM0290", "7", "HNF4A,NR2F1,RXRA"),
            tab_motif("21_5090724_A/T", "21:5090724", "T", "ENSM00000574815", "-1", "ENSPFM0484", "15", "RUNX3"),
            "21_5090724_A/T\t21:5090724\tT\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\t-\t-\t-\t-\t-".to_string(),
            "rs868224642\t21:46699255\tC\t-\tENSR21_C6M55\tRegulatoryFeature\tregulatory_region_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\tenhancer\t-\t-\t-\t-\t-".to_string(),
            tab_motif("rs868224642", "21:46699255", "C", "ENSM00000588223", "-1", "ENSPFM0427", "5", "NRF1"),
            tab_motif("rs868224642", "21:46699255", "C", "ENSM00000588224", "1", "ENSPFM0427", "8", "NRF1"),
            tab_motif("rs868224642", "21:46699255", "C", "ENSM00000588225", "-1", "ENSPFM0427", "11", "NRF1"),
            tab_motif("rs868224642", "21:46699255", "C", "ENSM00000588226", "1", "ENSPFM0427", "2", "NRF1"),
            "rs868224642\t21:46699255\tC\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\t-\t-\t-\t-\t-".to_string(),
        ]
    );
    let vcf: Vec<String> = rows("vcf")
        .iter()
        .map(|line| line.split('\t').nth(7).unwrap().to_string())
        .collect();
    let csq_motif = |allele: &str,
                     id: &str,
                     strand: &str,
                     matrix: &str,
                     pos: &str,
                     factors: &str| {
        format!("{allele}|TF_binding_site_variant|MODIFIER|||MotifFeature|{id}|||||||||||||{strand}||||{matrix}|{pos}|||{factors}")
    };
    assert_eq!(
        vcf,
        [
            format!(
                "CSQ={},{},T|intergenic_variant|MODIFIER|||||||||||||||||||||||||",
                csq_motif("T", "ENSM00000574814", "1", "ENSPFM0290", "7", "HNF4A&NR2F1&RXRA"),
                csq_motif("T", "ENSM00000574815", "-1", "ENSPFM0484", "15", "RUNX3")
            ),
            format!(
                "CSQ=C|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR21_C6M55|enhancer||||||||||||||||||||,{},{},{},{},C|intergenic_variant|MODIFIER|||||||||||||||||||||||||",
                csq_motif("C", "ENSM00000588223", "-1", "ENSPFM0427", "5", "NRF1"),
                csq_motif("C", "ENSM00000588224", "1", "ENSPFM0427", "8", "NRF1"),
                csq_motif("C", "ENSM00000588225", "-1", "ENSPFM0427", "11", "NRF1"),
                csq_motif("C", "ENSM00000588226", "1", "ENSPFM0427", "2", "NRF1")
            ),
        ]
    );
    let json: Vec<Value> = rows("json")
        .iter()
        .map(|line| {
            let mut record: Value = serde_json::from_str(line).unwrap();
            // The input line and the version key are vep-rs's own, outside the comparison.
            record.as_object_mut().unwrap().remove("input");
            record.as_object_mut().unwrap().remove("vep_rs_version");
            record
        })
        .collect();
    let motif = |id: &str, matrix: &str, pos: u32, strand: i32, factors: &[&str], allele: &str| {
        serde_json::json!({
            "consequence_terms": ["TF_binding_site_variant"], "impact": "MODIFIER",
            "motif_feature_id": id, "motif_name": matrix, "motif_pos": pos, "strand": strand,
            "transcription_factors": factors, "variant_allele": allele
        })
    };
    assert_eq!(
        json,
        [
            serde_json::json!({
                "allele_string": "A/T", "assembly_name": "GRCh38", "end": 5090724, "id": ".",
                "intergenic_consequences": [{"consequence_terms": ["intergenic_variant"], "impact": "MODIFIER", "variant_allele": "T"}],
                "most_severe_consequence": "TF_binding_site_variant",
                "motif_feature_consequences": [
                    motif("ENSM00000574814", "ENSPFM0290", 7, 1, &["HNF4A", "NR2F1", "RXRA"], "T"),
                    motif("ENSM00000574815", "ENSPFM0484", 15, -1, &["RUNX3"], "T"),
                ],
                "seq_region_name": "21", "start": 5090724, "strand": 1
            }),
            serde_json::json!({
                "allele_string": "G/C", "assembly_name": "GRCh38", "end": 46699255, "id": "rs868224642",
                "intergenic_consequences": [{"consequence_terms": ["intergenic_variant"], "impact": "MODIFIER", "variant_allele": "C"}],
                "most_severe_consequence": "TF_binding_site_variant",
                "motif_feature_consequences": [
                    motif("ENSM00000588223", "ENSPFM0427", 5, -1, &["NRF1"], "C"),
                    motif("ENSM00000588224", "ENSPFM0427", 8, 1, &["NRF1"], "C"),
                    motif("ENSM00000588225", "ENSPFM0427", 11, -1, &["NRF1"], "C"),
                    motif("ENSM00000588226", "ENSPFM0427", 2, 1, &["NRF1"], "C"),
                ],
                "regulatory_feature_consequences": [{"biotype": "enhancer", "consequence_terms": ["regulatory_region_variant"], "impact": "MODIFIER", "regulatory_feature_id": "ENSR21_C6M55", "variant_allele": "C"}],
                "seq_region_name": "21", "start": 46699255, "strand": 1
            }),
        ]
    );
}

/// The Parquet output holds the regulatory, motif and intergenic rows the
/// `--tab` output holds, in both shapes, read back through the adapter and
/// compared as a multiset (Parquet is sorted by variant key). A motif row's
/// factor list comes back as the tab cell; a factorless motif's empty factor
/// cell, a written value the tab format leaves empty rather than `-`, is NULL
/// in Parquet like every absent value and comes back as `-`. Runs only where
/// the DuckDB CLI the adapter needs is installed.
#[test]
fn regulatory_rows_round_trip_through_parquet() {
    let duckdb = std::process::Command::new("duckdb")
        .arg("--version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false);
    if !duckdb {
        eprintln!("skipping: duckdb CLI not installed");
        return;
    }
    let adapter =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../../scripts/adapters/parquet_to_vep_tab.py");
    let dir = tempfile::tempdir().unwrap();
    let cache = cache(dir.path());
    fs::write(
        dir.path().join("in.vcf"),
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t106\tforward\tC\tG\t.\t.\t.\n\
         21\t116\treverse\tA\tG\t.\t.\t.\n",
    )
    .unwrap();
    let tab = annotate(dir.path(), &cache, "tab", true);
    let mut tab_rows: Vec<String> = tab
        .lines()
        .filter(|l| !l.starts_with("##"))
        .map(|l| {
            if l.ends_with('\t') {
                format!("{l}-")
            } else {
                l.to_string()
            }
        })
        .collect();
    tab_rows.sort_unstable();
    assert!(
        tab_rows.iter().any(|row| row.contains("\tTF_A,TF_B")),
        "the forward motif row carries its factors:\n{tab}"
    );
    assert!(
        tab.lines()
            .any(|row| row.contains("ENSM_REVERSE") && row.ends_with('\t')),
        "the factorless motif row ends in the empty factor cell:\n{tab}"
    );
    for shape in ["flat", "nested"] {
        let parquet_dir = dir.path().join(format!("out-{shape}.parquet"));
        let mut command = common::isolated_vep(dir.path());
        command
            .args([
                "--json_cache",
                cache.to_str().unwrap(),
                "--assembly",
                "GRCh38",
                "--no_stats",
                "--no_fasta",
                "--force_overwrite",
                "--format",
                "vcf",
                "--regulatory",
                "--output_format",
                "parquet",
                "--parquet_shape",
                shape,
            ])
            .arg("--input_file")
            .arg(dir.path().join("in.vcf"))
            .arg("--output_file")
            .arg(&parquet_dir);
        let result = command.output().unwrap();
        assert!(
            result.status.success(),
            "parquet ({shape}): {}",
            String::from_utf8_lossy(&result.stderr)
        );
        let back = std::process::Command::new("python3")
            .arg(&adapter)
            .arg(&parquet_dir)
            .output()
            .unwrap();
        assert!(
            back.status.success(),
            "adapter ({shape}): {}",
            String::from_utf8_lossy(&back.stderr)
        );
        let round_trip = String::from_utf8(back.stdout).unwrap();
        let mut parquet_rows: Vec<String> = round_trip
            .lines()
            .filter(|l| !l.starts_with("##"))
            .map(str::to_string)
            .collect();
        parquet_rows.sort_unstable();
        assert_eq!(parquet_rows, tab_rows, "parquet ({shape}) round trip");
    }
}

/// Every structural class and the sequence-variant edges over the hand-built
/// fixture (enhancer 100-120, motifs 105-108 and 115-118), as Ensembl VEP
/// 116.2 writes them: `<DEL>` and `<CN0>` covering a feature are its ablation
/// (`deletion`, `copy_number_loss`: Utils/VariationEffect.pm:192-216), `<DUP>`,
/// `<CN2>`, `<DUP:TANDEM>` and `<CNV:TR>` its amplification (`copy_number_gain`,
/// :246-256), `<CNV>` and `<INV>` the overlap term alone (a VCF `<CNV>` has no
/// supporting variants), each with `OverlapBP` and `OverlapPC`
/// (OutputFactory.pm:2087-2095) and none of `BIOTYPE`, `STRAND` and the motif
/// fields (:2069-2102 writes no such key); `--overlaps` adds the two as columns
/// to the tab and VCF outputs, `-` and empty on the intergenic row. A sequence
/// deletion covering the features is their ablation with the motif fields and
/// no `MOTIF_POS` (its start lies before the matrix); a multi-allelic SNV writes
/// each allele's rows; an insertion overlaps a feature only with both flanks
/// inside it (InputBuffer.pm:324): 100-101 and 104-105 overlap the enhancer and
/// not the motif at 105-108, 99-100 and 120-121 overlap nothing.
#[test]
fn sv_and_indel_classes_over_a_regulatory_feature() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache(dir.path());
    let header = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
    let symbolic = [
        "<DEL>",
        "<DUP>",
        "<CNV>",
        "<INV>",
        "<CN0>",
        "<CN2>",
        "<CNV:TR>",
        "<DUP:TANDEM>",
    ];
    let mut input = header.to_string();
    for alt in symbolic {
        input.push_str(&format!("21\t94\t{alt}\tN\t{alt}\t.\t.\tEND=125\n"));
    }
    input.push_str("21\t104\tpartial\tN\t<DEL>\t.\t.\tEND=112\n");
    input.push_str("21\t99\tcover_del\tAAAAAAAAAAAAAAAAAAAAAA\tA\t.\t.\t.\n");
    input.push_str("21\t106\tmulti\tC\tG,T\t.\t.\t.\n");
    for (id, pos) in [
        ("ins_in", 100),
        ("ins_before", 99),
        ("ins_after", 120),
        ("ins_motif", 104),
    ] {
        input.push_str(&format!("21\t{pos}\t{id}\tA\tAT\t.\t.\t.\n"));
    }
    fs::write(dir.path().join("in.vcf"), input).unwrap();

    // (id, class term, regulatory terms, motif terms, motif impact)
    let classes = [
        (
            "<DEL>",
            "deletion",
            "regulatory_region_ablation,regulatory_region_variant",
            "TFBS_ablation,TF_binding_site_variant",
            "MODERATE",
        ),
        (
            "<DUP>",
            "duplication",
            "regulatory_region_amplification,regulatory_region_variant",
            "TFBS_amplification,TF_binding_site_variant",
            "MODIFIER",
        ),
        (
            "<CNV>",
            "copy_number_variation",
            "regulatory_region_variant",
            "TF_binding_site_variant",
            "MODIFIER",
        ),
        (
            "<INV>",
            "inversion",
            "regulatory_region_variant",
            "TF_binding_site_variant",
            "MODIFIER",
        ),
        (
            "<CN0>",
            "deletion",
            "regulatory_region_ablation,regulatory_region_variant",
            "TFBS_ablation,TF_binding_site_variant",
            "MODERATE",
        ),
        (
            "<CN2>",
            "duplication",
            "regulatory_region_amplification,regulatory_region_variant",
            "TFBS_amplification,TF_binding_site_variant",
            "MODIFIER",
        ),
        (
            "<CNV:TR>",
            "tandem_repeat",
            "regulatory_region_amplification,regulatory_region_variant",
            "TFBS_amplification,TF_binding_site_variant",
            "MODIFIER",
        ),
        (
            "<DUP:TANDEM>",
            "tandem_duplication",
            "regulatory_region_amplification,regulatory_region_variant",
            "TFBS_amplification,TF_binding_site_variant",
            "MODIFIER",
        ),
    ];
    let mut expected = Vec::new();
    for (id, class, regulatory, motif, motif_impact) in classes {
        let row = |feature: &str, kind: &str, terms: &str, extra: &str| {
            format!("{id}\t21:95-125\t{class}\t-\t{feature}\t{kind}\t{terms}\t-\t-\t-\t-\t-\t-\t{extra}")
        };
        expected.push(row(
            "ENSR_FIXTURE",
            "RegulatoryFeature",
            regulatory,
            "IMPACT=MODIFIER;OverlapBP=21;OverlapPC=100.00",
        ));
        for feature in ["ENSM_FORWARD", "ENSM_REVERSE"] {
            expected.push(row(
                feature,
                "MotifFeature",
                motif,
                &format!("IMPACT={motif_impact};OverlapBP=4;OverlapPC=100.00"),
            ));
        }
        expected.push(row("-", "-", "intergenic_variant", "IMPACT=MODIFIER"));
    }
    expected.extend([
        "partial\t21:105-112\tdeletion\t-\tENSR_FIXTURE\tRegulatoryFeature\tregulatory_region_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;OverlapBP=8;OverlapPC=38.10".to_string(),
        "partial\t21:105-112\tdeletion\t-\tENSM_FORWARD\tMotifFeature\tTFBS_ablation,TF_binding_site_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODERATE;OverlapBP=4;OverlapPC=100.00".to_string(),
        "partial\t21:105-112\tdeletion\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER".to_string(),
        "cover_del\t21:100-120\t-\t-\tENSR_FIXTURE\tRegulatoryFeature\tregulatory_region_ablation,regulatory_region_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;BIOTYPE=enhancer".to_string(),
        "cover_del\t21:100-120\t-\t-\tENSM_FORWARD\tMotifFeature\tTFBS_ablation,TF_binding_site_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODERATE;STRAND=1;MOTIF_NAME=ENSPFM_FORWARD;TRANSCRIPTION_FACTORS=TF_A,TF_B".to_string(),
        "cover_del\t21:100-120\t-\t-\tENSM_REVERSE\tMotifFeature\tTFBS_ablation,TF_binding_site_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODERATE;STRAND=-1;MOTIF_NAME=ENSPFM_REVERSE;TRANSCRIPTION_FACTORS=".to_string(),
        "cover_del\t21:100-120\t-\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER".to_string(),
        "multi\t21:106\tG\t-\tENSR_FIXTURE\tRegulatoryFeature\tregulatory_region_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;BIOTYPE=enhancer".to_string(),
        "multi\t21:106\tT\t-\tENSR_FIXTURE\tRegulatoryFeature\tregulatory_region_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;BIOTYPE=enhancer".to_string(),
        "multi\t21:106\tG\t-\tENSM_FORWARD\tMotifFeature\tTF_binding_site_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;MOTIF_NAME=ENSPFM_FORWARD;MOTIF_POS=2;TRANSCRIPTION_FACTORS=TF_A,TF_B".to_string(),
        "multi\t21:106\tT\t-\tENSM_FORWARD\tMotifFeature\tTF_binding_site_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;STRAND=1;MOTIF_NAME=ENSPFM_FORWARD;MOTIF_POS=2;TRANSCRIPTION_FACTORS=TF_A,TF_B".to_string(),
        "multi\t21:106\tG\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER".to_string(),
        "multi\t21:106\tT\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER".to_string(),
        "ins_in\t21:100-101\tT\t-\tENSR_FIXTURE\tRegulatoryFeature\tregulatory_region_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;BIOTYPE=enhancer".to_string(),
        "ins_in\t21:100-101\tT\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER".to_string(),
        "ins_before\t21:99-100\tT\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER".to_string(),
        "ins_after\t21:120-121\tT\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER".to_string(),
        "ins_motif\t21:104-105\tT\t-\tENSR_FIXTURE\tRegulatoryFeature\tregulatory_region_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER;BIOTYPE=enhancer".to_string(),
        "ins_motif\t21:104-105\tT\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tIMPACT=MODIFIER".to_string(),
    ]);
    let rows: Vec<String> = annotate(dir.path(), &cache, "default", true)
        .lines()
        .filter(|line| !line.starts_with('#'))
        .map(str::to_owned)
        .collect();
    assert_eq!(rows, expected);

    // --overlaps: the two values become columns of the tab and VCF outputs, after
    // FLAGS and before BIOTYPE, `-` and empty on the intergenic row.
    let tab = annotate_with_args(dir.path(), &cache, "tab", true, &["--overlaps"]);
    let header = tab
        .lines()
        .find(|l| l.starts_with("#Uploaded_variation"))
        .unwrap();
    assert_eq!(
        header,
        "#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\tConsequence\tcDNA_position\tCDS_position\tProtein_position\tAmino_acids\tCodons\tExisting_variation\tIMPACT\tDISTANCE\tSTRAND\tFLAGS\tOverlapBP\tOverlapPC\tBIOTYPE\tMOTIF_NAME\tMOTIF_POS\tHIGH_INF_POS\tMOTIF_SCORE_CHANGE\tTRANSCRIPTION_FACTORS"
    );
    let partial: Vec<&str> = tab.lines().filter(|l| l.starts_with("partial\t")).collect();
    assert_eq!(
        partial,
        [
            "partial\t21:105-112\tdeletion\t-\tENSR_FIXTURE\tRegulatoryFeature\tregulatory_region_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t8\t38.10\t-\t-\t-\t-\t-\t-",
            "partial\t21:105-112\tdeletion\t-\tENSM_FORWARD\tMotifFeature\tTFBS_ablation,TF_binding_site_variant\t-\t-\t-\t-\t-\t-\tMODERATE\t-\t-\t-\t4\t100.00\t-\t-\t-\t-\t-\t-",
            "partial\t21:105-112\tdeletion\t-\t-\t-\tintergenic_variant\t-\t-\t-\t-\t-\t-\tMODIFIER\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-\t-",
        ]
    );
    let vcf = annotate_with_args(dir.path(), &cache, "vcf", true, &["--overlaps"]);
    let format = vcf
        .lines()
        .find(|l| l.starts_with("##INFO=<ID=CSQ"))
        .and_then(|l| l.split("Format: ").nth(1))
        .unwrap();
    assert_eq!(
        format,
        "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|OverlapBP|OverlapPC|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS\">"
    );
    let csq = |id: &str| -> Vec<String> {
        vcf.lines()
            .find(|l| l.split('\t').nth(2) == Some(id))
            .and_then(|l| l.split('\t').nth(7))
            .and_then(|info| info.split(';').find_map(|item| item.strip_prefix("CSQ=")))
            .unwrap()
            .split(',')
            .map(str::to_owned)
            .collect()
    };
    assert_eq!(
        csq("partial"),
        [
            "deletion|regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR_FIXTURE|||||||||||||||||8|38.10|||||",
            "deletion|TFBS_ablation&TF_binding_site_variant|MODERATE|||MotifFeature|ENSM_FORWARD|||||||||||||||||4|100.00|||||",
            "deletion|intergenic_variant|MODIFIER|||||||||||||||||||||||||||",
        ]
    );
    assert_eq!(
        csq("cover_del"),
        [
            "-|regulatory_region_ablation&regulatory_region_variant|MODIFIER|||RegulatoryFeature|ENSR_FIXTURE|enhancer||||||||||||||||||||||",
            "-|TFBS_ablation&TF_binding_site_variant|MODERATE|||MotifFeature|ENSM_FORWARD|||||||||||||1||||||ENSPFM_FORWARD||||TF_A&TF_B",
            "-|TFBS_ablation&TF_binding_site_variant|MODERATE|||MotifFeature|ENSM_REVERSE|||||||||||||-1||||||ENSPFM_REVERSE||||",
            "-|intergenic_variant|MODIFIER|||||||||||||||||||||||||||",
        ]
    );
}

/// The pick family through the command line, on a record inside two
/// transcripts of one gene and the fixture's features shifted onto it, and on a
/// record inside the same features shifted between the genes: `--pick` keeps
/// the one worst allele of the record (the intron row; the motif row where no
/// transcript overlaps, the intergenic row then dropped); `--per_gene` and
/// `--pick_allele_gene` keep one transcript row per gene and every regulatory,
/// motif and intergenic row (OutputFactory.pm:837-845, non-transcript alleles
/// pass through); `--flag_pick` and `--flag_pick_allele` flag the one pick;
/// `--flag_pick_allele_gene` flags the per-gene pick and every pass-through row,
/// the intergenic row included (:623-627); `--most_severe` and `--summary` read
/// the regulatory terms with the transcript terms.
#[test]
fn pick_family_flags_through_the_cli_keep_regulatory_rows() {
    let dir = tempfile::tempdir().unwrap();
    let cache = dir.path().join("cache");
    fs::create_dir_all(cache.join("transcripts/21")).unwrap();
    fs::create_dir_all(cache.join("regulatory/21")).unwrap();
    // Two protein-coding transcripts of ENSG00000166351 at 13609777-13645823 and
    // three lncRNA transcripts of another gene at 13843132-13897298.
    fs::copy(
        Path::new(env!("CARGO_MANIFEST_DIR")).join(
            "../../tests/golden/116/GRCh38/json_cache/transcripts/21/13000001-14000000.json.gz",
        ),
        cache.join("transcripts/21/13000001-14000000.json.gz"),
    )
    .unwrap();
    let fixture: Vec<Value> =
        serde_json::from_slice(include_bytes!("fixtures/regulatory/1-1000000.json")).unwrap();
    let shifted: Vec<Value> = [(13_619_894u64, "_IN_GENE"), (13_699_894, "_BETWEEN")]
        .iter()
        .flat_map(|(offset, suffix)| {
            fixture.iter().map(move |feature| {
                let mut feature = feature.clone();
                for key in ["start", "end"] {
                    let value = feature[key].as_u64().unwrap() + offset;
                    feature[key] = Value::from(value);
                }
                let id = format!("{}{suffix}", feature["stable_id"].as_str().unwrap());
                feature["stable_id"] = Value::from(id);
                feature
            })
        })
        .collect();
    fs::write(
        cache.join("regulatory/21/13000001-14000000.json"),
        serde_json::to_vec(&shifted).unwrap(),
    )
    .unwrap();
    fs::write(
        cache.join("info.json"),
        r#"{"assembly":"GRCh38","regulatory":true}"#,
    )
    .unwrap();
    fs::write(
        dir.path().join("in.vcf"),
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
         21\t13620000\tin_gene\tC\tG\t.\t.\t.\n\
         21\t13700000\tbetween\tC\tG\t.\t.\t.\n",
    )
    .unwrap();
    // (record, Feature, Feature_type, Consequence, PICK)
    let rows = |flags: &[&str]| -> Vec<(String, String, String, String, String)> {
        annotate_with_args(dir.path(), &cache, "default", true, flags)
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(|line| {
                let cols: Vec<&str> = line.split('\t').collect();
                let pick = cols[13]
                    .split(';')
                    .find_map(|item| item.strip_prefix("PICK="))
                    .unwrap_or("-");
                (
                    cols[0].to_string(),
                    cols[4].to_string(),
                    cols[5].to_string(),
                    cols[6].to_string(),
                    pick.to_string(),
                )
            })
            .collect()
    };
    let row = |record: &str, feature: &str, kind: &str, terms: &str, pick: &str| {
        (
            record.to_string(),
            feature.to_string(),
            kind.to_string(),
            terms.to_string(),
            pick.to_string(),
        )
    };
    let in_gene_t1 = |pick| {
        row(
            "in_gene",
            "ENST00000299443",
            "Transcript",
            "intron_variant",
            pick,
        )
    };
    let in_gene_t2 = |pick| {
        row(
            "in_gene",
            "ENST00000620442",
            "Transcript",
            "intron_variant",
            pick,
        )
    };
    let in_gene_rf = |pick| {
        row(
            "in_gene",
            "ENSR_FIXTURE_IN_GENE",
            "RegulatoryFeature",
            "regulatory_region_variant",
            pick,
        )
    };
    let in_gene_mf = |pick| {
        row(
            "in_gene",
            "ENSM_FORWARD_IN_GENE",
            "MotifFeature",
            "TF_binding_site_variant",
            pick,
        )
    };
    let between_rf = |pick| {
        row(
            "between",
            "ENSR_FIXTURE_BETWEEN",
            "RegulatoryFeature",
            "regulatory_region_variant",
            pick,
        )
    };
    let between_mf = |pick| {
        row(
            "between",
            "ENSM_FORWARD_BETWEEN",
            "MotifFeature",
            "TF_binding_site_variant",
            pick,
        )
    };
    let between_ig = |pick| row("between", "-", "-", "intergenic_variant", pick);

    assert_eq!(
        rows(&[]),
        [
            in_gene_t1("-"),
            in_gene_t2("-"),
            in_gene_rf("-"),
            in_gene_mf("-"),
            between_rf("-"),
            between_mf("-"),
            between_ig("-")
        ]
    );
    assert_eq!(rows(&["--pick"]), [in_gene_t1("-"), between_mf("-")]);
    for flag in ["--per_gene", "--pick_allele_gene"] {
        assert_eq!(
            rows(&[flag]),
            [
                in_gene_t1("-"),
                in_gene_rf("-"),
                in_gene_mf("-"),
                between_rf("-"),
                between_mf("-"),
                between_ig("-")
            ],
            "{flag}"
        );
    }
    for flag in ["--flag_pick", "--flag_pick_allele"] {
        assert_eq!(
            rows(&[flag]),
            [
                in_gene_t1("1"),
                in_gene_t2("-"),
                in_gene_rf("-"),
                in_gene_mf("-"),
                between_rf("-"),
                between_mf("1"),
                between_ig("-")
            ],
            "{flag}"
        );
    }
    assert_eq!(
        rows(&["--flag_pick_allele_gene"]),
        [
            in_gene_t1("1"),
            in_gene_t2("-"),
            in_gene_rf("1"),
            in_gene_mf("1"),
            between_rf("1"),
            between_mf("1"),
            between_ig("1")
        ]
    );
    assert_eq!(
        rows(&["--most_severe"]),
        [
            row("in_gene", "-", "-", "intron_variant", "-"),
            row("between", "-", "-", "TF_binding_site_variant", "-")
        ]
    );
    assert_eq!(
        rows(&["--summary"]),
        [
            row(
                "in_gene",
                "-",
                "-",
                "intron_variant,TF_binding_site_variant,regulatory_region_variant",
                "-"
            ),
            row(
                "between",
                "-",
                "-",
                "TF_binding_site_variant,regulatory_region_variant,intergenic_variant",
                "-"
            ),
        ]
    );
}

/// A regulatory shard the run reads must be a JSON array of regulatory and
/// motif features: a truncated shard, a feature type other than the two (a
/// transcript among them), an element without the feature fields and a motif
/// without its binding matrix each fail the run with exit 1 and the shard's
/// path in the error, writing no data row; an empty array is a chromosome
/// with no features, so its records write their intergenic rows; a
/// gzip-compressed shard loads as a plain one does.
#[test]
fn malformed_regulatory_shards_fail_the_run_by_name() {
    let dir = tempfile::tempdir().unwrap();
    fs::write(
        dir.path().join("in.vcf"),
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t106\tr\tC\tG\t.\t.\t.\n",
    )
    .unwrap();
    let run = |shard: &[u8], name: &str| -> (std::process::ExitStatus, String, Vec<String>) {
        let cache = dir.path().join("cache");
        let _ = fs::remove_dir_all(&cache);
        fs::create_dir_all(cache.join("transcripts/21")).unwrap();
        fs::create_dir_all(cache.join("regulatory/21")).unwrap();
        fs::write(cache.join("transcripts/21/1-1000000.json"), "[]").unwrap();
        fs::write(cache.join("regulatory/21").join(name), shard).unwrap();
        fs::write(
            cache.join("info.json"),
            r#"{"assembly":"GRCh38","regulatory":true}"#,
        )
        .unwrap();
        let out = dir.path().join("out.txt");
        let _ = fs::remove_file(&out);
        let result = common::isolated_vep(dir.path())
            .args(["--json_cache", cache.to_str().unwrap()])
            .args([
                "--assembly",
                "GRCh38",
                "--no_fasta",
                "--no_stats",
                "--force_overwrite",
            ])
            .arg("--regulatory")
            .arg("--input_file")
            .arg(dir.path().join("in.vcf"))
            .arg("--output_file")
            .arg(&out)
            .output()
            .unwrap();
        let rows = fs::read_to_string(&out)
            .unwrap_or_default()
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(str::to_owned)
            .collect();
        (
            result.status,
            String::from_utf8_lossy(&result.stderr).into_owned(),
            rows,
        )
    };
    for (shard, phrase) in [
        (&b"["[..], "EOF while parsing"),
        (
            br#"[{"stable_id":"X","start":100,"end":120,"strand":0,"feature_type":"Enhancer"}]"#,
            "unknown variant `Enhancer`, expected `RegulatoryFeature` or `MotifFeature`",
        ),
        (
            br#"[{"stable_id":"T1","start":100,"end":120,"strand":1,"feature_type":"Transcript"}]"#,
            "unknown variant `Transcript`",
        ),
        (br#"[{"foo":1}]"#, "missing field `stable_id`"),
        (
            br#"[{"stable_id":"ENSM_X","start":100,"end":120,"strand":1,"feature_type":"MotifFeature"}]"#,
            "motif feature ENSM_X has no binding matrix",
        ),
    ] {
        let (status, stderr, rows) = run(shard, "1-1000000.json");
        assert_eq!(status.code(), Some(1), "{stderr}");
        assert!(
            stderr.contains("ERROR: [cache] parsing")
                && stderr.contains("regulatory/21/1-1000000.json")
                && stderr.contains(phrase),
            "{stderr}"
        );
        assert!(rows.is_empty(), "{rows:?}");
    }
    let (status, stderr, rows) = run(b"[]", "1-1000000.json");
    assert!(status.success(), "{stderr}");
    assert_eq!(rows.len(), 1);
    assert!(
        rows[0].contains("\t-\t-\t-\tintergenic_variant\t"),
        "{}",
        rows[0]
    );
    let mut encoder = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::default());
    std::io::Write::write_all(
        &mut encoder,
        include_bytes!("fixtures/regulatory/1-1000000.json"),
    )
    .unwrap();
    let (status, stderr, rows) = run(&encoder.finish().unwrap(), "1-1000000.json.gz");
    assert!(status.success(), "{stderr}");
    assert_eq!(rows.len(), 3, "{rows:?}");
    assert!(
        rows[0].contains("\tENSR_FIXTURE\tRegulatoryFeature\t"),
        "{}",
        rows[0]
    );
}

#[test]
fn regulatory_and_motif_rows_use_dash_for_absent_gene() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache(dir.path());
    fs::write(
        dir.path().join("in.vcf"),
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t106\tforward\tC\tG\t.\t.\t.\n",
    )
    .unwrap();
    for format in ["default", "tab"] {
        let text = annotate(dir.path(), &cache, format, true);
        let features: Vec<_> = text
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(|line| line.split('\t').collect::<Vec<_>>())
            .filter(|columns| matches!(columns[5], "RegulatoryFeature" | "MotifFeature"))
            .map(|columns| (columns[4], columns[3]))
            .collect();
        assert_eq!(
            features,
            [("ENSR_FIXTURE", "-"), ("ENSM_FORWARD", "-")],
            "{format}: {text}"
        );
    }
}

fn assert_regulatory_only_rows(extra_args: &[&str], regulatory: bool, intergenic: bool) {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache(dir.path());
    fs::write(
        dir.path().join("in.vcf"),
        "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n21\t110\tregulatory_only\tA\tG\t.\t.\t.\n",
    )
    .unwrap();
    for format in ["default", "tab"] {
        let text = annotate_with_args(dir.path(), &cache, format, true, extra_args);
        let rows: Vec<_> = text
            .lines()
            .filter(|line| !line.starts_with('#'))
            .map(|line| line.split('\t').collect::<Vec<_>>())
            .collect();
        let mut expected = Vec::new();
        if regulatory {
            expected.push((
                "ENSR_FIXTURE",
                "RegulatoryFeature",
                "regulatory_region_variant",
            ));
        }
        if intergenic {
            expected.push(("-", "-", "intergenic_variant"));
        }
        assert_eq!(
            rows.iter()
                .map(|row| (row[4], row[5], row[6]))
                .collect::<Vec<_>>(),
            expected,
            "{format} {extra_args:?}: {text}"
        );
        for row in &rows {
            assert_eq!(row[3], "-");
            if row[5] == "RegulatoryFeature" {
                if format == "default" {
                    assert!(row[13].split(';').any(|item| item == "BIOTYPE=enhancer"));
                } else {
                    let header = text
                        .lines()
                        .find(|line| line.starts_with("#Uploaded_variation"))
                        .unwrap();
                    let biotype = header
                        .split('\t')
                        .position(|column| column == "BIOTYPE")
                        .unwrap();
                    assert_eq!(row[biotype], "enhancer");
                }
            }
        }
    }
    let text = annotate_with_args(dir.path(), &cache, "json", true, extra_args);
    let records: Vec<Value> = text
        .lines()
        .map(|line| serde_json::from_str(line).unwrap())
        .collect();
    assert_eq!(records.len(), 1);
    let record = &records[0];
    assert!(record.get("transcript_consequences").is_none());
    assert!(record.get("motif_feature_consequences").is_none());
    assert_eq!(
        record["most_severe_consequence"],
        if regulatory {
            "regulatory_region_variant"
        } else {
            "?"
        }
    );
    if regulatory {
        assert_eq!(
            record["regulatory_feature_consequences"],
            serde_json::json!([{
                "biotype": "enhancer",
                "consequence_terms": ["regulatory_region_variant"],
                "impact": "MODIFIER",
                "regulatory_feature_id": "ENSR_FIXTURE",
                "variant_allele": "G"
            }])
        );
    } else {
        assert!(record.get("regulatory_feature_consequences").is_none());
    }
    if intergenic {
        assert_eq!(
            record["intergenic_consequences"],
            serde_json::json!([{
                "consequence_terms": ["intergenic_variant"],
                "impact": "MODIFIER",
                "variant_allele": "G"
            }])
        );
    } else {
        assert!(record.get("intergenic_consequences").is_none());
    }
    let text = annotate_with_args(dir.path(), &cache, "vcf", true, extra_args);
    let rows: Vec<_> = text.lines().filter(|line| !line.starts_with('#')).collect();
    assert_eq!(rows.len(), 1);
    let info = rows[0].split('\t').nth(7).unwrap();
    if regulatory {
        let consequences: Vec<_> = info.strip_prefix("CSQ=").unwrap().split(',').collect();
        assert_eq!(consequences.len(), 1 + usize::from(intergenic));
        let regulatory_row: Vec<_> = consequences[0].split('|').collect();
        assert_eq!(
            &regulatory_row[..8],
            &[
                "G",
                "regulatory_region_variant",
                "MODIFIER",
                "",
                "",
                "RegulatoryFeature",
                "ENSR_FIXTURE",
                "enhancer"
            ]
        );
        if intergenic {
            let intergenic_row: Vec<_> = consequences[1].split('|').collect();
            assert_eq!(
                &intergenic_row[..8],
                &["G", "intergenic_variant", "MODIFIER", "", "", "", "", ""]
            );
        }
    } else {
        assert_eq!(info, ".");
    }
}

#[test]
fn regulatory_only_retains_intergenic_row() {
    assert_regulatory_only_rows(&[], true, true);
}

#[test]
fn regulatory_only_pick_retains_regulatory_row() {
    assert_regulatory_only_rows(&["--pick"], true, false);
}

#[test]
fn regulatory_only_coding_only_removes_all_rows() {
    assert_regulatory_only_rows(&["--coding_only"], false, false);
}

#[test]
fn regulatory_only_canonical_retains_intergenic_and_regulatory_rows() {
    assert_regulatory_only_rows(&["--canonical"], true, true);
}

#[test]
fn regulatory_only_biotype_retains_intergenic_and_regulatory_rows() {
    assert_regulatory_only_rows(&["--biotype"], true, true);
}

#[test]
fn regulatory_only_no_intergenic_removes_all_rows() {
    assert_regulatory_only_rows(&["--no_intergenic"], false, false);
}

/// A corrupt regulatory shard fails the embedding annotator's construction, as a
/// corrupt transcript shard does, instead of annotating every variant of the
/// chromosome without its regulatory rows and reporting the failure as an event.
#[test]
fn embedding_annotator_refuses_a_corrupt_regulatory_shard_at_construction() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache(dir.path());
    let config = vep_cli::Config::from_argv([
        "vep",
        "--json_cache",
        cache.to_str().unwrap(),
        "--regulatory",
        "--no_fasta",
    ])
    .unwrap();
    let error = match vep_cli::Annotator::from_config(&config) {
        Err(error) => format!("{error:#}"),
        Ok(_) => panic!("a cache with a corrupt regulatory shard built an annotator"),
    };
    assert!(
        error.contains("parsing") && error.contains("regulatory/22/1-1000000.json"),
        "{error}"
    );
}

#[test]
fn embedding_batches_load_regulatory_features_and_apply_summary_filters() {
    let dir = tempfile::tempdir().unwrap();
    let cache = cache(dir.path());
    // The embedding annotator loads every chromosome when it is built.
    fs::remove_dir_all(cache.join("regulatory/22")).unwrap();
    let config = vep_cli::Config::from_argv([
        "vep",
        "--json_cache",
        cache.to_str().unwrap(),
        "--regulatory",
        "--no_fasta",
    ])
    .unwrap();
    let annotator = vep_cli::Annotator::from_config(&config).unwrap();
    let mut variants = vec![
        vep_cli::InputVariant::new("21".into(), 106, 106, b"C".to_vec(), b"G".to_vec()),
        vep_cli::InputVariant::new("21".into(), 115, 115, b"A".to_vec(), b"G".to_vec()),
    ];
    annotator.annotate_batch(&mut variants);
    assert_eq!(
        variants[0].most_severe_consequence,
        Some(vep_cli::Consequence::TfBindingSiteVariant)
    );
    assert!(variants[0].intergenic_overlap);
    assert!(variants[1]
        .transcript_consequences
        .iter()
        .any(|r| r.motif.as_deref().is_some_and(|m| m.pos == Some(4))));
    let mut summary = variants.clone();
    vep_cli::apply_to_batch(
        &mut summary,
        &vep_cli::FilterConfig {
            summary: true,
            ..Default::default()
        },
    );
    assert_eq!(
        summary[0].summary_terms.as_deref(),
        Some(
            &[
                vep_cli::Consequence::TfBindingSiteVariant,
                vep_cli::Consequence::RegulatoryRegionVariant,
                vep_cli::Consequence::IntergenicVariant,
            ][..]
        )
    );
    let mut coding = variants.clone();
    vep_cli::apply_to_batch(
        &mut coding,
        &vep_cli::FilterConfig {
            row_filters: vep_cli::RowFilters {
                coding_only: true,
                ..Default::default()
            },
            ..Default::default()
        },
    );
    assert!(coding
        .iter()
        .all(|v| v.transcript_consequences.is_empty() && v.most_severe_consequence.is_none()));
    vep_cli::apply_to_batch(
        &mut variants,
        &vep_cli::FilterConfig {
            row_filters: vep_cli::RowFilters {
                no_intergenic: true,
                ..Default::default()
            },
            ..Default::default()
        },
    );
    assert!(variants
        .iter()
        .all(|v| v.transcript_consequences.is_empty() && v.most_severe_consequence.is_none()));
}
